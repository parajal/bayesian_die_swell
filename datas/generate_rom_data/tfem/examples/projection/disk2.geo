ox = -0.5;
oy = -0.5;
lx = 1.0;
ly = 1.0;
dx_box = 0.1;
nobj = 1;
xp[1] = 0.0;
yp[1] = 0.0;
rp[1] = 0.1;
dx_part = 0.017;

Mesh.Algorithm = 6; // Frontal 
Mesh.ElementOrder = 2; // Second-order elements

Include "particles_in_a_box_2D.igo";
