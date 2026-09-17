L = 15;
H = 2;
a = 1;
//M0:
hcyl   = 0.08;
hwall  = 0.4;
hLwall = 1.2;
hLcen  = 0.08;
//M1:
//hcyl   = 0.04;
//hwall  = 0.2;
//hLwall = 0.6;
//hLcen  = 0.04;
//M2:
//hcyl   = 0.02;
//hwall  = 0.1;
//hLwall = 0.3;
//hLcen  = 0.02;
//M3:
//hcyl   = 0.01;
//hwall  = 0.05;
//hLwall = 0.15;
//hLcen  = 0.01;
//M4:
//hcyl   = 0.005;
//hwall  = 0.025;
//hLwall = 0.075;
//hLcen  = 0.005;

Include "confined_cylinder.igo";

Mesh.Algorithm = 6; // Frontal
Mesh.ElementOrder = 2; // Second-order elements

