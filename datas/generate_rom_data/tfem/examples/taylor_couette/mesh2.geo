// Geometric parameters

x[1] =   0.5; 
y[1] =   0.0;
x[2] =   1.0;
y[2] =   0.0;

theta =  2.0; // angle of slice (NOTE: should be smaller than pi)


// meshing size

n_th = 40; // number of element in angle theta


// Create points and line
 
Point(1) = { x[1], y[1],  0.0 };
Point(2) = { x[2], y[2],  0.0 };

Line(1) = {1, 2};


// Extrude the line

out[] = Extrude {{0,0,1}, {0,0,0}, theta} {
     Line{1}; Layers{n_th};
};


// Define physical points, lines and surface

Physical Point(1) = 1;
Physical Point(2) = 2;
Physical Point(3) = 3;
Physical Point(4) = 4;

Physical Curve("cut_1",   1) = 1;
Physical Curve("cut_2",   2) = 2;
Physical Curve("wall_in",  3) = 3;
Physical Curve("wall_out", 4) = 4;

Physical Surface("cavity", 1) = 5;


// Transfinite for structured quad mesh

Transfinite Line {1, 2} = 41 Using Progression 1;
Transfinite Surface {5} = {1, 2, 3, 4};
Recombine Surface {5};


// Meshing

Mesh 2;
SetOrder 2;
RenumberMeshNodes;
RenumberMeshElements;
Mesh.MshFileVersion = 2.2;

